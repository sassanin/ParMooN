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
// @(#)DirectSolver.h
//
// Purpose:     solve equation system by direct solver
//
// Author:      Gunar Matthies (06.09.05)
//
// History:     start of implementation 06.09.05 (Gunar Matthies)
//
// =======================================================================

#include <MainUtilities.h>
#include <DirectSolver.h>
#include <Database.h>
#include "stdlib.h"
#include <LinAlg.h>

extern "C"
{
  #include "umfpack.h"
}
void UMFPACK_return(int ret)
{

  if (TDatabase::ParamDB->SC_VERBOSE>=2)
  {
      OutPut("solved " << ret << endl);
  }
  if (ret == 0) 
  {
     if (TDatabase::ParamDB->SC_VERBOSE>=2)
      OutPut("solved sucessfully (" << ret << ")" << endl);
  }
  else if (ret ==   1) 
  {
      OutPut("solved failed (" << ret << "): Matrix singular" << endl);
  }
  else if (ret ==   2)
  {
      OutPut("solved failed (" << ret << "): det(A) != 0 but < eps" << endl);
  }
  else if (ret ==   3)
  {
      OutPut("solved failed (" << ret << "): det(A) != 0 but > inf" << endl);
  }
  else if (ret ==  -1)
  {
      OutPut("solved failed (" << ret << "): Not enough memory!" << endl);
  }
  else if (ret ==  -3)
  {
      OutPut("solved failed (" << ret << "): Used Numeric object is invalided!" << endl);
  }
  else if (ret ==  -4)
  {
      OutPut("solved failed (" << ret << "): Used Symbolic object is invalided!" << endl);
  }
  else if (ret ==  -5)
  {
      OutPut("solved failed (" << ret << "): Argument missing!" << endl);
  }
  else if (ret ==  -6)
  {
      OutPut("solved failed (" << ret << "): Number of rows and columns must be greater 0!" << endl);
  }
  else if (ret ==  -8)
  {
      OutPut("solved failed (" << ret << "): Invalidid matrix!" << endl);
  }
  else if (ret == -11)
  {
      OutPut("solved failed (" << ret << "): Different pattern" << endl);
  }
  else if (ret == -13)
  {
      OutPut("solved failed (" << ret << "): Invalidid system!" << endl);
  }
  else if (ret == -15)
  {
      OutPut("solved failed (" << ret << "): Invalidid permutation!" << endl);
  }
  else if (ret == -17)
  {
      OutPut("solved failed (" << ret << "): Error due to I/O!!!" << endl);
  }
  else if (ret == -911)
  {
      OutPut("solved failed (" << ret << "): Internal error!!!" << endl);
  }

  if (ret != 0) exit(4711);
}

/*******************************************************************/
/*        SCALAR PROBLEMS                                          */
/*******************************************************************/
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn;
  int *Row, *KCol;
  double *Values;
  void *Symbolic, *Numeric;

  N_Eqn = matrix->GetN_Columns();
  Row = matrix->GetRowPtr();
  KCol = matrix->GetKCol();
  Values = matrix->GetEntries();

  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    OutPut("umfpack: reordering of the columns will be performed"<<endl);
    OutPut("umfpack: no back ordering implemented !!!"<<endl);

    for(i=0;i<N_Eqn;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Values[j];
            KCol[j] = KCol[k]; Values[j] = Values[k];
            KCol[k] = l;       Values[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
  }
 
  t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values,
    &Symbolic, NULL, NULL);
  t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_symbolic " << ret << endl);
    exit(4711);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic,
    &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  t3 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_numeric " << ret << endl);
    exit(4711);
  }

  ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values,
    sol, rhs, Numeric, NULL, NULL);
  umfpack_di_free_numeric(&Numeric);
  t4 = GetTime();
  if (ret!=0)
  {
    OutPut("error in umfpack_di_solve " << ret << endl);
    exit(4711);
  }
//  OutPut("umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}

/*******************************************************************/
/*        SCALAR PROBLEMS                                          */
/*******************************************************************/
void DirectSolver(TSquareMatrix2D *matrix, double *rhs, double *sol)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn;
  int *Row, *KCol;
  double *Values;
  void *Symbolic, *Numeric;

  N_Eqn = matrix->GetN_Columns();
  Row = matrix->GetRowPtr();
  KCol = matrix->GetKCol();
  Values = matrix->GetEntries();

  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    OutPut("umfpack: reordering of the columns will be performed"<<endl);
    OutPut("umfpack: no back ordering implemented !!!"<<endl);

    for(i=0;i<N_Eqn;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Values[j];
            KCol[j] = KCol[k]; Values[j] = Values[k];
            KCol[k] = l;       Values[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
  }

  t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values,
    &Symbolic, NULL, NULL);
  t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_symbolic " << ret << endl);
    exit(4711);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic,
    &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  t3 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_numeric " << ret << endl);
    exit(4711);
  }

  ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values,
    sol, rhs, Numeric, NULL, NULL);

  umfpack_di_free_numeric(&Numeric);

  t4 = GetTime();
  if (ret!=0)
  {
    OutPut("error in umfpack_di_solve " << ret << endl);
    exit(4711);
  }
 OutPut("umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}


// rb_flag = 0 ==> allocation and LU-decomposition forward/backward.
// rb_flag = 1 ==> only forward/backward.
// rb_flag = 2 ==> forward/backward and free up memory
// rb_flag = 3 ==> allocation, LU-decomposition, forward/backward, free up memory
// rb_flag = 4 ==> only free up memory
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn, N_Entries;
  int *Row_orig, *KCol_orig;
  double *Values_orig;
//   void *Symbolic, *Numeric;

//   static double *Values;
//   static int *KCol, *Row;
//   static void *Symbolic, *Numeric;  
//     
  double *null = (double *) NULL;

  
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
    return;
  }  
  
  if (TDatabase::ParamDB->SC_VERBOSE>=3)
  {
   OutPut("rb_flag: " << rb_flag << endl);
  }
  
  
 if (rb_flag==0 || rb_flag==3)
 {  
  N_Eqn = matrix->GetN_Columns();
  Row_orig = matrix->GetRowPtr();
  KCol_orig = matrix->GetKCol();
  Values_orig = matrix->GetEntries();
  
  N_Entries = Row_orig[N_Eqn];     
  KCol = new int[N_Entries];
  Row = new int[N_Eqn+1];
  Values = new double[N_Entries];
  
  
  memcpy(Values, Values_orig, N_Entries*SizeOfDouble);
  memcpy(KCol, KCol_orig, N_Entries*SizeOfInt);
  memcpy(Row, Row_orig, (N_Eqn+1)*SizeOfInt);
   
   
  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    OutPut("umfpack: reordering of the columns will be performed"<<endl);
    OutPut("umfpack: no back ordering implemented !!!"<<endl);

    for(i=0;i<N_Eqn;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Values[j];
            KCol[j] = KCol[k]; Values[j] = Values[k];
            KCol[k] = l;       Values[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
  }
 
  t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values, &Symbolic, NULL, NULL);
  t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_symbolic " << ret << endl);
    exit(4711);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic, &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  t3 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_numeric " << ret << endl);
    exit(4711);
  }
 } // if (rb_flag==0 || rb_flag==3)
 
  
  ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values, sol, rhs, Numeric, NULL, NULL);
//   umfpack_di_free_numeric(&Numeric);
  t4 = GetTime();
  if (ret!=0)
  {
    OutPut("error in umfpack_di_solve " << ret << endl);
    exit(4711);
  }
  
  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
  }
  
//  OutPut("umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}




/*******************************************************************/
/*        SCALAR PROBLEMS with multiple rhs                        */
/*******************************************************************/
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn;
  int *Row, *KCol;
  double *Values, *Sol, *Rhs;
  void *Symbolic, *Numeric;

  N_Eqn = matrix->GetN_Columns();
  Row = matrix->GetRowPtr();
  KCol = matrix->GetKCol();
  Values = matrix->GetEntries();

  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    OutPut("umfpack: reordering of the columns will be performed"<<endl);
    OutPut("umfpack: no back ordering implemented !!!"<<endl);

    for(i=0;i<N_Eqn;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Values[j];
            KCol[j] = KCol[k]; Values[j] = Values[k];
            KCol[k] = l;       Values[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
  }

  //t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values, &Symbolic, NULL, NULL);
  //t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_symbolic " << ret << endl);
    //exit(4711);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic, &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  //t3 = GetTime();
  // error occured
  if (ret!=0)
   {
    OutPut("error in umfpack_di_numeric " << ret << endl);
    //exit(4711);
   }

 for(i=N_Rhs_Disp; i<N_Rhs; i++)
  {
   //Sol = sol[i];
   Sol = sol+i*N_Eqn;
   Rhs = rhs+i*N_Eqn;

   ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values, Sol, Rhs, Numeric, NULL, NULL);

   //t4 = GetTime();
   if (ret!=0)
    {
     OutPut("error in umfpack_di_solve " << ret << endl);
    //exit(4711);
    }
  } // for(i=0; i<N_Rhs; i++)

  umfpack_di_free_numeric(&Numeric);

//  OutPut("umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}



// rb_flag = 0 ==> allocation and LU-decomposition forward/backward.
// rb_flag = 1 ==> only forward/backward.
// rb_flag = 2 ==> forward/backward and free up memory
// rb_flag = 3 ==> allocation, LU-decomposition, forward/backward, free up memory
// rb_flag = 4 ==> only free up memory
/*******************************************************************/
/*        SCALAR PROBLEMS with multiple rhs                        */
/*******************************************************************/
void DirectSolver(TSquareMatrix *matrix, double *rhs, double *sol, int N_Rhs, int N_Rhs_Disp, double *&Values,
                   int *&KCol, int *&Row, void *&Symbolic, void *&Numeric, int rb_flag)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  int N_Eqn, N_Entries;
  int *Row_orig, *KCol_orig;
  double *Values_orig, *Sol, *Rhs;
//   void *Symbolic, *Numeric;
 
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
    return;
  }    
 
 
  if (TDatabase::ParamDB->SC_VERBOSE>=3)
  {
   OutPut("rb_flag: " << rb_flag << endl);
  }
  
  N_Eqn = matrix->GetN_Columns(); 
 if (rb_flag==0 || rb_flag==3)
 {   
  Row_orig = matrix->GetRowPtr();
  KCol_orig = matrix->GetKCol();
  Values_orig = matrix->GetEntries();

  N_Entries = Row_orig[N_Eqn];     
  KCol = new int[N_Entries];
  Row = new int[N_Eqn+1];
  Values = new double[N_Entries];
    
  memcpy(Values, Values_orig, N_Entries*SizeOfDouble);
  memcpy(KCol, KCol_orig, N_Entries*SizeOfInt);
  memcpy(Row, Row_orig, (N_Eqn+1)*SizeOfInt);
  
  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
    // sort matrix
    OutPut("umfpack: reordering of the columns will be performed"<<endl);
    OutPut("umfpack: no back ordering implemented !!!"<<endl);

    for(i=0;i<N_Eqn;i++)
    {
      begin=Row[i];
      end=Row[i+1];
      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];     value = Values[j];
            KCol[j] = KCol[k]; Values[j] = Values[k];
            KCol[k] = l;       Values[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i
  }

  //t1 = GetTime();
  ret = umfpack_di_symbolic(N_Eqn, N_Eqn, Row, KCol, Values, &Symbolic, NULL, NULL);
  //t2 = GetTime();
  // error occured
  if (ret!=0)
  {
    OutPut("error in umfpack_di_symbolic " << ret << endl);
    //exit(4711);
  }

  ret = umfpack_di_numeric(Row, KCol, Values, Symbolic, &Numeric, NULL, NULL);
  umfpack_di_free_symbolic(&Symbolic);
  //t3 = GetTime();
  // error occured
  if (ret!=0)
   {
    OutPut("error in umfpack_di_numeric " << ret << endl);
    //exit(4711);
   }
 }//if (rb_flag==0 || rb_flag==3)
 
 
 for(i=N_Rhs_Disp; i<N_Rhs; i++)
  {
   //Sol = sol[i];
   Sol = sol+i*N_Eqn;
   Rhs = rhs+i*N_Eqn;

   ret = umfpack_di_solve(UMFPACK_At, Row, KCol, Values, Sol, Rhs, Numeric, NULL, NULL);

   //t4 = GetTime();
   if (ret!=0)
    {
     OutPut("error in umfpack_di_solve " << ret << endl);
    //exit(4711);
    }
  } // for(i=0; i<N_Rhs; i++)

//   umfpack_di_free_numeric(&Numeric);

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Values;
    delete [] KCol;
    delete [] Row;
  }

//  OutPut("umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}



void DirectSolverLong(TSquareMatrix *matrix, double *rhs, double *sol)
{
  double t1, t2, t3, t4;
  int ret, i, j, k, l, begin, end;
  double value;
  long N_Eqn;
  long *Row, *KCol;
  int *row, *kcol;
  double *Values;
  void *Symbolic, *Numeric;

  N_Eqn = matrix->GetN_Columns();
  row = matrix->GetRowPtr();
  kcol = matrix->GetKCol();
  Values = matrix->GetEntries();

  
  Row = new long[N_Eqn+1];
  for(i=0;i<=N_Eqn;i++)
    Row[i] = row[i];

  end = Row[N_Eqn];
  KCol = new long[end];
  for(i=0;i<end;i++)
    KCol[i] = kcol[i];  
  
  // check ordering of the matrix
  if (matrix->GetColOrder() != 1)
  {
      // sort matrix
      OutPut("umfpack: reordering of the columns will be performed"<<endl);
      OutPut("umfpack: no back ordering implemented !!!"<<endl);

      for(i=0;i<N_Eqn;i++)
      {
          begin=Row[i];
          end=Row[i+1];
          for(j=begin;j<end;j++)
          {
              for(k=j+1;k<end;k++)
              {
                  if(KCol[j] > KCol[k])
                  {
                      l = KCol[j];     value = Values[j];
                      KCol[j] = KCol[k]; Values[j] = Values[k];
                      KCol[k] = l;       Values[k] = value;
                  } // endif
              } // endfor k
          } // endfor j
      } // endfor i
  }

  t1 = GetTime();
  ret = umfpack_dl_symbolic(N_Eqn, N_Eqn, Row, KCol, Values,
                          &Symbolic, NULL, NULL);
  t2 = GetTime();
//   OutPut("symbolic: " << ret << " " << t2-t1 << endl);
 
  ret = umfpack_dl_numeric(Row, KCol, Values, Symbolic,
                          &Numeric, NULL, NULL);
  umfpack_dl_free_symbolic(&Symbolic);
  t3 = GetTime();
//   OutPut("numeric: " << ret << " "  << t3-t2 << endl);
 
  ret = umfpack_dl_solve(UMFPACK_At, Row, KCol, Values,
                       sol, rhs, Numeric, NULL, NULL);
  umfpack_dl_free_numeric(&Numeric);
  delete [] KCol;
  delete [] Row;
  t4 = GetTime();
//   OutPut("solve: " << ret << " " << t4-t3 << endl);
//   OutPut("long umfpack: " << ret << " " << t4-t1 << " sec." << endl);
}


// Solver for PDAE2D
//
// rb_flag = 0 ==> allocieren und LU-Zerlegung
// rb_flag = 1 ==> only forward/backward.
// rb_flag = 2 ==> forward/backward and free up memory
// rb_flag = 3 ==> allocieren, LU-Zerl, Freigabe
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
double *rhs1, double *rhs2, double *sol1, double *sol2, int rb_flag)
{
  int *KColA, *RowPtrA;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *sol, *rhs;
  int N_, N_U, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  N_U = sqmatrixA11->GetN_Rows();
  N_ = 2*N_U;

  // copy sol and rhs piece by piece
  sol = new double[N_];
  rhs = new double[N_];
  for (i=0;i<N_U;i++)
  {
    sol[i] = sol1[i];
    rhs[i] = rhs1[i];
  }
  for (i=N_U;i<N_;i++)
  {
    sol[i] = sol2[i-N_U];
    rhs[i] = rhs2[i-N_U];
  }

  OutPut("rb_flag: " << rb_flag << endl);
  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    N_Active = sqmatrixA11->GetActiveBound();

    KColA = sqmatrixA11->GetKCol();
    RowPtrA = sqmatrixA11->GetRowPtr();

    EntriesA11 = sqmatrixA11->GetEntries();
    EntriesA12 = sqmatrixA12->GetEntries();
    EntriesA21 = sqmatrixA21->GetEntries();
    EntriesA22 = sqmatrixA22->GetEntries();

    N_Entries = 4*RowPtrA[N_U];
    Entries = new double[N_Entries];
    KCol = new int[N_Entries];

    RowPtr = new int[N_+1];
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

      RowPtr[N_U+i+1] = pos;
    }

    // sort matrix
    for(i=0;i<N_;i++)
    {
      begin=RowPtr[i];
      end=RowPtr[i+1];

      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];      value = Entries[j];
            KCol[j] = KCol[k]; Entries[j] = Entries[k];
            KCol[k] = l;       Entries[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i

    /*
    for(i=0;i<N_;i++)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        cout << i << " " << KCol[j] << " " << Entries[j] << endl;
    }
    */

//     t2 = GetTime();
    ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//     OutPut("symbolic: " << ret << endl);
//     t3 = GetTime();

    ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//     OutPut("numeric: " << ret << endl);
//     t4 = GetTime();
    umfpack_di_free_symbolic(&Symbolic);
  }

  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
//   OutPut("solve: " << ret << endl);
  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }

  // copy sol and rhs piece by piece
  for (i=0;i<N_U;i++)
    sol1[i] = sol[i];
  for (i=N_U;i<N_;i++)
    sol2[i-N_U] = sol[i];

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;

  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}


// rb_flag = 0 ==> allocieren und LU-Zerlegung
// rb_flag = 1 ==> nur vorw./rueckw.
// rb_flag = 2 ==> speicher wieder freigeben
// rb_flag = 3 ==> allocieren, LU-Zerl, Freigabe
void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
double *rhs, double *sol, int rb_flag)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;
  int verbose = TDatabase::ParamDB->SC_VERBOSE;
  double sum = 0;

  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
    return;
  }

  if (TDatabase::ParamDB->SC_VERBOSE>=3)
  {
      OutPut("rb_flag: " << rb_flag << endl);
  }

  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    N_U = sqmatrixA11->GetN_Rows();
    N_P = matrixB1->GetN_Rows();
    N_ = 2*N_U + N_P;
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
    Entries = new double[N_Entries];
    KCol = new int[N_Entries];

    RowPtr = new int[N_+1];
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

    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      // pressure constant
      begin = RowPtr[2*N_U];
      end = RowPtr[2*N_U+1];
      for(j=begin+1;j<end;j++)
        Entries[j] = 0;
      Entries[begin] = 1;
      KCol[begin] = 2*N_U;
      rhs[2*N_U] = 0;
    }

    // sort matrix
    for(i=0;i<N_;i++)
    {
      begin=RowPtr[i];
      end=RowPtr[i+1];

      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];      value = Entries[j];
            KCol[j] = KCol[k]; Entries[j] = Entries[k];
            KCol[k] = l;       Entries[k] = value;
          }                      // endif
        }                        // endfor k
      }                          // endfor j
    }                            // endfor i

    /*
    for(i=0;i<N_;i++)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        cout << i << " " << KCol[j] << " " << Entries[j] << endl;
    }
    */

    t2 = GetTime();
    ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
    if (ret!=0)
    {
        OutPut("WARNING: symbolic: " << ret << endl);
    }
    t3 = GetTime();

    ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
    if (ret!=0)
    {
        OutPut("WARNING: numeric: " << ret << endl);
    }
    t4 = GetTime();
    umfpack_di_free_symbolic(&Symbolic);
  }

  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if (ret!=0)
  {
      OutPut("WARNING: solve: " << ret << endl);
  }

  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }

  if (verbose>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }

  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}

void DirectSolver(TSquareMatrix2D *sqmatrixA,
TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
TMatrix2D *matrixC,
double *rhs, double *sol)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  int *KColC, *RowPtrC;
  double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T, *EntriesC;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();

  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();

  KColC = matrixC->GetKCol();
  RowPtrC = matrixC->GetRowPtr();

  EntriesA = sqmatrixA->GetEntries();
  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();
  EntriesC = matrixC->GetEntries();

  N_Entries = 2*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U] + RowPtrC[N_P];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
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

    begin = RowPtrC[i];
    end = RowPtrC[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesC[j];
      KCol[pos] = KColC[j] + 2*N_U;
      pos++;
    }

    RowPtr[2*N_U+i+1] = pos;
  }

  /*
    // check residual
    {
      double *res, val, sum;
      res = new double[N_];
      sum = 0;
      for(i=0;i<N_;i++)
      {
        val = rhs[i];
        begin = RowPtr[i];
        end = RowPtr[i+1];
        for(j=begin;j<end;j++)
          val -= Entries[j]*sol[KCol[j]];
        res[i] = val;
        sum += val*val;
        if(fabs(val)>1e-12)
        {
          cout << setw(5) << i << setw(5) << i-2*N_U << setw(25) << res[i] << setw(25) << rhs[i] << endl;
        }
      }
      cout << "defect: " << sqrt(sum) << endl;
    }
  */

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[2*N_U];
    end = RowPtr[2*N_U+1];
    k = begin;
    for(j=begin+1;j<end;j++)
    {
      Entries[j] = 0;
      // check whether col (2*N_U) is already in this row
      if(KCol[j] == 2*N_U) k = j;
    }
    Entries[k] = 1;
    KCol[k] = 2*N_U;
    rhs[2*N_U] = 0;
  }

  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    cout << "row: " << i << " " << RowPtr[i] << ".." << RowPtr[i+1] << endl;
    for(j=RowPtr[i];j<RowPtr[i+1]-1;j++)
      if(KCol[j]>=KCol[j+1]) cout << "Error: " << j << endl;

    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//   OutPut("symbolic: " << ret << endl);
//   t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//   OutPut("numeric: " << ret << endl);
//   t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries, sol, rhs, Numeric, null, null);
//   OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;
}


//****************************************************************************/
//
// for NSTYPE == 1
//
//****************************************************************************/

void DirectSolver(TSquareMatrix2D *sqmatrixA,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
double *rhs, double *sol)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  double *EntriesA, *EntriesB1, *EntriesB2;
  int N_, N_U, N_P, N_B, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5, sum;

  t1 = GetTime();
  // get information from the matrices
  // size
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();
  // pointer to the index arrays
  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  // entries
  EntriesA = sqmatrixA->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();

  // allocate arrays for structure of combined matrix
  // total number of entries
  N_Entries = 2*RowPtrA[N_U] + 4*RowPtrB[N_P];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];
  RowPtr = new int[N_+1];
  RowPtr[0] = 0;
  N_B = RowPtrB[N_P];

  pos = 0;
  // fill combined matrix
  for(i=0;i<N_U;i++)
  {
    // first velocity component
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
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

  /*  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      // pressure constant
      begin = RowPtr[2*N_U];
      end = RowPtr[2*N_U+1];
      for(j=begin+1;j<end;j++)
        Entries[j] = 0;
      Entries[begin] = 1;
      KCol[begin] = 2*N_U;
      rhs[2*N_U] = 0;
    }
  */

  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//   OutPut("symbolic: " << ret << endl);
//   t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//   OutPut("numeric: " << ret << endl);
//   t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries, sol, rhs, Numeric, null, null);
//   OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;

  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}

void DirectSolver(TSquareMatrix2D *sqmatrixA,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
double *rhs, double *sol, int rb_flag)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  double *EntriesA, *EntriesB1, *EntriesB2;
  int N_, N_U, N_P, N_B, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5, sum;
  
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
      OutPut("rb_flag: " << rb_flag << endl);
    return;
  }

  if (TDatabase::ParamDB->SC_VERBOSE>=3)
  {
      OutPut("rb_flag: " << rb_flag << endl);
  }

  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    // get information from the matrices
    // size
    N_U = sqmatrixA->GetN_Rows();
    N_P = matrixB1->GetN_Rows();
    N_ = 2*N_U + N_P;
    N_Active = sqmatrixA->GetActiveBound();
    // pointer to the index arrays
    KColA = sqmatrixA->GetKCol();
    RowPtrA = sqmatrixA->GetRowPtr();

    KColB = matrixB1->GetKCol();
    RowPtrB = matrixB1->GetRowPtr();

    // entries
    EntriesA = sqmatrixA->GetEntries();

    EntriesB1 = matrixB1->GetEntries();
    EntriesB2 = matrixB2->GetEntries();

    // allocate arrays for structure of combined matrix
    // total number of entries
    N_Entries = 2*RowPtrA[N_U] + 4*RowPtrB[N_P];
    Entries = new double[N_Entries];
    KCol = new int[N_Entries];
    RowPtr = new int[N_+1];
    RowPtr[0] = 0;
    N_B = RowPtrB[N_P];

    pos = 0;
    // fill combined matrix
    for(i=0;i<N_U;i++)
    {
      // first velocity component
      begin = RowPtrA[i];
      end = RowPtrA[i+1];
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

    /*  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      // pressure constant
      begin = RowPtr[2*N_U];
      end = RowPtr[2*N_U+1];
      for(j=begin+1;j<end;j++)
        Entries[j] = 0;
      Entries[begin] = 1;
      KCol[begin] = 2*N_U;
      rhs[2*N_U] = 0;
    }
    */

    // sort matrix
    for(i=0;i<N_;i++)
    {
      begin=RowPtr[i];
      end=RowPtr[i+1];

      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];      value = Entries[j];
            KCol[j] = KCol[k]; Entries[j] = Entries[k];
            KCol[k] = l;       Entries[k] = value;
          }                        // endif
        }                          // endfor k
      }                            // endfor j
    }                              // endfor i

    /*
    for(i=0;i<N_;i++)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        cout << i << " " << KCol[j] << " " << Entries[j] << endl;
    }
    */

//     t2 = GetTime();

    ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//     OutPut("symbolic: " << ret << endl);
//     t3 = GetTime();
    ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//     OutPut("numeric: " << ret << endl);
//     t4 = GetTime();
    umfpack_di_free_symbolic(&Symbolic);
  }

  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if (TDatabase::ParamDB->SC_VERBOSE>=2)
  {
      OutPut("solve: " << ret << endl);
  }
  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }
  
//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;

  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}

//****************************************************************************/
//
// for NSTYPE == 2
//
//****************************************************************************/

void DirectSolver(TSquareMatrix2D *sqmatrixA,
TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
double *rhs, double *sol)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
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
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
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

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[2*N_U];
    end = RowPtr[2*N_U+1];
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 2*N_U;
    rhs[2*N_U] = 0;
  }

  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */

//   t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
//   OutPut("symbolic: " << ret << endl);
//   t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
//   OutPut("numeric: " << ret << endl);
//   t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
//   OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
//   t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;
}

void DirectSolver(TSquareMatrix2D *sqmatrixA,
TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
double *rhs, double *sol, int rb_flag)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA, *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  if (TDatabase::ParamDB->SC_VERBOSE>=3)
  {
      OutPut("rb_flag: " << rb_flag << endl);
  }
  if (rb_flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;

    return;
  }

  if (rb_flag==0 || rb_flag==3)
  {
    t1 = GetTime();
    N_U = sqmatrixA->GetN_Rows();
    N_P = matrixB1->GetN_Rows();
    N_ = 2*N_U + N_P;
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
    Entries = new double[N_Entries];
    KCol = new int[N_Entries];

    RowPtr = new int[N_+1];
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

    /*
    if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      // pressure constant
      begin = RowPtr[2*N_U];
      end = RowPtr[2*N_U+1];
      for(j=begin+1;j<end;j++)
        Entries[j] = 0;
      Entries[begin] = 1;
      KCol[begin] = 2*N_U;
      rhs[2*N_U] = 0;
    }
    */

    // sort matrix
    for(i=0;i<N_;i++)
    {
      begin=RowPtr[i];
      end=RowPtr[i+1];

      for(j=begin;j<end;j++)
      {
        for(k=j+1;k<end;k++)
        {
          if(KCol[j] > KCol[k])
          {
            l = KCol[j];      value = Entries[j];
            KCol[j] = KCol[k]; Entries[j] = Entries[k];
            KCol[k] = l;       Entries[k] = value;
          }                        // endif
        }                          // endfor k
      }                            // endfor j
    }                              // endfor i

    /*
    for(i=0;i<N_;i++)
    {
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        cout << i << " " << KCol[j] << " " << Entries[j] << endl;
    }
    */

    t2 = GetTime();

    ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
    UMFPACK_return(ret);

    t3 = GetTime();
    ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
    UMFPACK_return(ret);

    t4 = GetTime();
    umfpack_di_free_symbolic(&Symbolic);
  }

  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  UMFPACK_return(ret);

  t5 = GetTime();

  if (rb_flag==2 || rb_flag==3)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }
  
  if (TDatabase::ParamDB->SC_VERBOSE>=2)
  {
    cout << "UMFPACK time:";
    cout << "  data prep : " << t2-t1 << " ";
    cout << "  symbolic: " << t3-t2 << " ";
    cout << "  numeric: " << t4-t3 << " ";
    cout << "  solve: " << t5-t4 << endl;
    cout << "UMFPACK total time: " << t5-t1 << endl;
  }
}

//****************************************************************************/
//
// for NSTYPE == 4
//
//****************************************************************************/

void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
double *rhs, double *sol)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = sqmatrixA11->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
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
  
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
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

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[2*N_U];
    end = RowPtr[2*N_U+1];
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 2*N_U;
    rhs[2*N_U] = 0;
  }

  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i
/*
    double sum1 = 0.0, sum2 = 0.0;
    for(i=0;i<N_P;i++)
  {
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      sum1 += EntriesB1[j];
      sum2 += EntriesB2[j];
    }
  }
  cout<<"B1 : "<<sum1<<"\n"<<"B2 : "<<sum2<<"\n";
  
    sum1= 0.0; 
    sum2 = 0.0;
    
      for(i=0;i<N_U;i++)
  {
    begin = RowPtrBT[i];
    end = RowPtrBT[i+1];
    for(j=begin;j<end;j++)
    {
      sum1 += EntriesB1T[j];
      sum2 += EntriesB2T[j];
    }
  }
  cout<<"B1T : "<<sum1<<"\n"<<"B2T : "<<sum2<<"\n";
  
    double sum11= 0.0, sum12=0.0, sum21=0.0, sum22=0.0; 
    
      for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      sum11 += EntriesA11[j];
      sum12 += EntriesA12[j];
      sum21 += EntriesA21[j];
      sum22 += EntriesA22[j];
    }
  }
  cout<<"A11 : "<<sum11<<"\nA12 : "<<sum12<<"\nA21 : "<<sum21<<"\nA22 : "<<sum22<<"\n";
  
  
 cout<<"sol u: "<<Ddot(2*N_U,sol, sol)<<"\n";
 cout<<"sol p: "<<Ddot(N_P,sol+2*N_U, sol+2*N_U)<<"\n";
   cout<<"rhs u: "<<Ddot(2*N_U,rhs, rhs)<<"\n";
 cout<<"rhs p: "<<Ddot(N_P,rhs+2*N_U, rhs+2*N_U)<<"\n";*/

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  if(ret!=0)
  OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)  
  OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)    
  OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();
  
//   cout<<"\n"<<"After solve\n";
//   
//  cout<<"sol u: "<<Ddot(2*N_U,sol, sol)<<"\n";
//  cout<<"sol p: "<<Ddot(N_P,sol+2*N_U, sol+2*N_U)<<"\n";
//  cout<<"rhs u: "<<Ddot(2*N_U,rhs, rhs)<<"\n";
//  cout<<"rhs p: "<<Ddot(N_P,rhs+2*N_U, rhs+2*N_U)<<"\n";
//  exit(0);
//   
  
  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

  
  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }
  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}


//****************************************************************************/
//
// for NSTYPE == 14
//
//****************************************************************************/

void DirectSolver(TSquareMatrix2D *sqmatrixA11, TSquareMatrix2D *sqmatrixA12,
TSquareMatrix2D *sqmatrixA21, TSquareMatrix2D *sqmatrixA22,
TSquareMatrix2D *sqmatrixC,
TMatrix2D *matrixB1T, TMatrix2D *matrixB2T,
TMatrix2D *matrixB1,  TMatrix2D *matrixB2,
double *rhs, double *sol)
{
  int *KColA, *RowPtrA, *KColC, *RowPtrC;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22, *EntriesC;
  double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  // get information from the matrices
  // size
  N_U = sqmatrixA11->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 2*N_U + N_P;
  N_Active = sqmatrixA11->GetActiveBound();
  // pointer to the index arrays
  KColA = sqmatrixA11->GetKCol();
  RowPtrA = sqmatrixA11->GetRowPtr();

  KColC = sqmatrixC->GetKCol();
  RowPtrC = sqmatrixC->GetRowPtr();

  KColB = matrixB1->GetKCol();
  RowPtrB = matrixB1->GetRowPtr();

  KColBT = matrixB1T->GetKCol();
  RowPtrBT = matrixB1T->GetRowPtr();
  // entries
  EntriesA11 = sqmatrixA11->GetEntries();
  EntriesA12 = sqmatrixA12->GetEntries();
  EntriesA21 = sqmatrixA21->GetEntries();
  EntriesA22 = sqmatrixA22->GetEntries();
  EntriesC = sqmatrixC->GetEntries();

  EntriesB1 = matrixB1->GetEntries();
  EntriesB2 = matrixB2->GetEntries();
  EntriesB1T = matrixB1T->GetEntries();
  EntriesB2T = matrixB2T->GetEntries();

  // allocate arrays for structure of combined matrix
  // total number of entries
  N_Entries = 4*RowPtrA[N_U] + RowPtrC[N_P] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];
  RowPtr = new int[N_+1];
  RowPtr[0] = 0;

  pos = 0;
  // fill combined matrix
  for(i=0;i<N_U;i++)
  {
    // first velocity component
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
    }
    // B1T
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
  // second velocity component
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
    }
    // B2T
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
    // C
    begin = RowPtrC[i];
    end = RowPtrC[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesC[j];
      KCol[pos] = KColC[j]+2*N_U;
      pos++;
    }
    RowPtr[2*N_U+i+1] = pos;
  }

  /*  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
    {
      // pressure constant
      begin = RowPtr[2*N_U];
      end = RowPtr[2*N_U+1];
      for(j=begin+1;j<end;j++)
        Entries[j] = 0;
      Entries[begin] = 1;
      KCol[begin] = 2*N_U;
      rhs[2*N_U] = 0;
    }
  */

  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete  []RowPtr;

  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }
  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}

void DirectSolver(TSquareMatrix2D *sqmatrixA, TSquareMatrix2D *sqmatrixC,
                  TMatrix2D *matrixBT, TMatrix2D *matrixB,
                  double *rhs, double *sol)
{
  int *KColA, *RowPtrA, *KColC, *RowPtrC;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA, *EntriesC, *EntriesB, *EntriesBT;
  int N_DOF, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  // get information from the matrices
  // size
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB->GetN_Rows();
  N_DOF = N_U + N_P;
  N_Active = sqmatrixA->GetActiveBound();
  
  // pointer to the index arrays
  KColA = sqmatrixA->GetKCol();
  RowPtrA = sqmatrixA->GetRowPtr();

  KColC = sqmatrixC->GetKCol();
  RowPtrC = sqmatrixC->GetRowPtr();

  KColB = matrixB->GetKCol();
  RowPtrB = matrixB->GetRowPtr();

  KColBT = matrixBT->GetKCol();
  RowPtrBT = matrixBT->GetRowPtr();
  // entries
  EntriesA = sqmatrixA->GetEntries();
  EntriesC = sqmatrixC->GetEntries();

  EntriesB = matrixB->GetEntries();
  EntriesBT = matrixBT->GetEntries();
  
  // allocate arrays for structure of combined matrix
  // total number of entries
  N_Entries = 4*RowPtrA[N_U] + RowPtrC[N_P] + RowPtrB[N_P] + RowPtrBT[N_U];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];
  RowPtr = new int[N_DOF+1];
  RowPtr[0] = 0;

  pos = 0;
  // fill combined matrix
  for(i=0;i<N_U;i++)
  {
    // first velocity component
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      // A
      Entries[pos] = EntriesA[j];
      KCol[pos] = KColA[j];
      pos++;
    }
    // BT
    //if(i<N_Active)
    {
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesBT[j];
        KCol[pos] = KColBT[j]+N_U;
        pos++;
      }
    }
    RowPtr[i+1] = pos;
  }
  // pressure
  for(i=0;i<N_P;i++)
  {
    // B
    begin = RowPtrB[i];
    end = RowPtrB[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesB[j];
      KCol[pos] = KColB[j];
      pos++;
    }
    // C
    begin = RowPtrC[i];
    end = RowPtrC[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesC[j];
      KCol[pos] = KColC[j]+N_U;
      pos++;
    }
    RowPtr[N_U+i+1] = pos;
  }

  // sort matrix
  for(i=0;i<N_DOF;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_DOF;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << "A(" << i+1 << "," << KCol[j]+1 << ") = " << Entries[j] << ";\n";
  }
  for(i=0;i<N_DOF;i++)
  {
    OutPut("rhs(" << i+1 << ") = " << rhs[i] << endl);
  }
  */
  
  t2 = GetTime();

  ret = umfpack_di_symbolic(N_DOF, N_DOF, RowPtr, KCol, Entries, &Symbolic, 
                            null, null);
  if(ret!=0)
    OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, 
                           null);
  if(ret!=0)
    OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)
    OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }
}



// for Conformation Stress Tensor
void DirectSolver(TSquareMatrix2D *sqmatrixS11, TSquareMatrix2D *sqmatrixS12, 
		  TSquareMatrix2D *sqmatrixS21, TSquareMatrix2D *sqmatrixS22, TSquareMatrix2D *sqmatrixS23, 
		  TSquareMatrix2D *sqmatrixS32, TSquareMatrix2D *sqmatrixS33,
                  double *rhs, double *sol)
{
  int *KColS, *RowPtrS;
  double *EntriesS11, *EntriesS12, *EntriesS21, *EntriesS22, *EntriesS23, *EntriesS32, *EntriesS33;
  int N_, N_S, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_S = sqmatrixS11->GetN_Rows();
  N_ = 3*N_S;
  N_Active = sqmatrixS11->GetActiveBound();

  KColS = sqmatrixS11->GetKCol();
  RowPtrS = sqmatrixS11->GetRowPtr();

  EntriesS11 = sqmatrixS11->GetEntries();
  EntriesS12 = sqmatrixS12->GetEntries();
  EntriesS21 = sqmatrixS21->GetEntries();
  EntriesS22 = sqmatrixS22->GetEntries();
  EntriesS23 = sqmatrixS23->GetEntries();
  EntriesS32 = sqmatrixS32->GetEntries();
  EntriesS33 = sqmatrixS33->GetEntries();
  
  N_Entries = 7*RowPtrS[N_S];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
  RowPtr[0] = 0;

  pos = 0;

  for(i=0;i<N_S;i++)
  {
    begin = RowPtrS[i];
    end = RowPtrS[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesS11[j];
      KCol[pos] = KColS[j];
      pos++;

      Entries[pos] = (i<N_Active)?EntriesS12[j]:0;
      KCol[pos] = KColS[j]+N_S;
      pos++;
    }
    RowPtr[i+1] = pos;
  }

  for(i=0;i<N_S;i++)
  {
    begin = RowPtrS[i];
    end = RowPtrS[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active)?EntriesS21[j]:0;
      KCol[pos] = KColS[j];
      pos++;

      Entries[pos] = EntriesS22[j];
      KCol[pos] = KColS[j]+N_S;
      pos++;
      
       Entries[pos] = (i<N_Active)?EntriesS23[j]:0;
      KCol[pos] = KColS[j]+(2*N_S);
      pos++;
      
    }
   RowPtr[N_S+i+1] = pos;
  }

    for(i=0;i<N_S;i++)
  {
    begin = RowPtrS[i];
    end = RowPtrS[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active)?EntriesS32[j]:0;
      KCol[pos] = KColS[j]+N_S;
      pos++;
      
      Entries[pos] = EntriesS33[j];
      KCol[pos] = KColS[j]+(2*N_S);
      pos++;
      
    }
   RowPtr[(2*N_S)+i+1] = pos;
  }



  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */
  
  
  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  if(ret!=0)
  OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)  
  OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)    
  OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();
// exit(0);
  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }
}

// for Deformation Tensor in DEVSS
void DirectSolver(TSquareMatrix2D *sqmatrixS11, TSquareMatrix2D *sqmatrixS22, 
		   TSquareMatrix2D *sqmatrixS33, double *rhs, double *sol)
{
  int *KColS, *RowPtrS;
  double *EntriesS11, *EntriesS22, *EntriesS33;
  int N_, N_S, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_S = sqmatrixS11->GetN_Rows();
  N_ = 3*N_S;
  N_Active = sqmatrixS11->GetActiveBound();

  KColS = sqmatrixS11->GetKCol();
  RowPtrS = sqmatrixS11->GetRowPtr();

  EntriesS11 = sqmatrixS11->GetEntries();
  EntriesS22 = sqmatrixS22->GetEntries();
  EntriesS33 = sqmatrixS33->GetEntries();
  
  N_Entries = 3*RowPtrS[N_S];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
  RowPtr[0] = 0;

  pos = 0;

  for(i=0;i<N_S;i++)
  {
    begin = RowPtrS[i];
    end = RowPtrS[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesS11[j];
      KCol[pos] = KColS[j];
      pos++;
    }
    RowPtr[i+1] = pos;
  }

  for(i=0;i<N_S;i++)
  {
    begin = RowPtrS[i];
    end = RowPtrS[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesS22[j];
      KCol[pos] = KColS[j]+N_S;
      pos++;
           
    }
   RowPtr[N_S+i+1] = pos;
  }

    for(i=0;i<N_S;i++)
  {
    begin = RowPtrS[i];
    end = RowPtrS[i+1];
    for(j=begin;j<end;j++)
    {
     
      Entries[pos] = EntriesS33[j];
      KCol[pos] = KColS[j]+(2*N_S);
      pos++;
      
    }
   RowPtr[(2*N_S)+i+1] = pos;
  }



  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */
  
  
  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  if(ret!=0)
  OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)  
  OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)    
  OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();
// exit(0);
  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }
}


// NSE(Type-4)-CST-2D (Monolithic) DEVSS
void DirectSolver(TSquareMatrix2D *SqmatrixA11, TSquareMatrix2D *SqmatrixA12, TSquareMatrix2D *SqmatrixA21, TSquareMatrix2D *SqmatrixA22,
		  TSquareMatrix2D *SqmatrixG11, TSquareMatrix2D *SqmatrixG12, TSquareMatrix2D *SqmatrixG21 , TSquareMatrix2D *SqmatrixG22, TSquareMatrix2D *SqmatrixG23 , TSquareMatrix2D *SqmatrixG32 , TSquareMatrix2D *SqmatrixG33,
                  TSquareMatrix2D *SqmatrixH11, TSquareMatrix2D *SqmatrixH22, TSquareMatrix2D *SqmatrixH33,
		  TMatrix2D *MatrixB1,  TMatrix2D *MatrixB2, TMatrix2D *MatrixB1T, TMatrix2D *MatrixB2T, 
                  TMatrix2D *MatrixC11, TMatrix2D *MatrixC12, TMatrix2D *MatrixC22, TMatrix2D *MatrixC23,
		  TMatrix2D *MatrixD11, TMatrix2D *MatrixD12, TMatrix2D *MatrixD21, TMatrix2D *MatrixD22, TMatrix2D *MatrixD31, TMatrix2D *MatrixD32,
	          TMatrix2D *MatrixE11, TMatrix2D *MatrixE12, TMatrix2D *MatrixE22, TMatrix2D *MatrixE23,
	          TMatrix2D *MatrixJ11, TMatrix2D *MatrixJ21, TMatrix2D *MatrixJ22, TMatrix2D *MatrixJ32,
                  double *rhs, double *sol)
{
  int *KColA, *RowPtrA, *KColB, *RowPtrB, *KColBT, *RowPtrBT, *KColG, *RowPtrG, *KColH, *RowPtrH;
  int *KColC, *RowPtrC, *KColD, *RowPtrD, *KColE, *RowPtrE, *KColJ, *RowPtrJ;
  
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG22, *EntriesG23, *EntriesG32, *EntriesG33;
  double *EntriesH11, *EntriesH22, *EntriesH33; 
  double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  double *EntriesC11, *EntriesC12, *EntriesC22, *EntriesC23;
  double *EntriesD11, *EntriesD12, *EntriesD21, *EntriesD22, *EntriesD31, *EntriesD32;
  double *EntriesE11, *EntriesE12, *EntriesE22, *EntriesE23;
  double *EntriesJ11, *EntriesJ21, *EntriesJ22, *EntriesJ32;
  
  int N_Tot, N_U, N_P, N_S, N_D, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active_U, N_Active_S, N_Active_D;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = SqmatrixA11->GetN_Rows();
  N_S = SqmatrixG11->GetN_Rows();
  N_D = SqmatrixH11->GetN_Rows();
  N_P = MatrixB1->GetN_Rows();
  N_Tot = 2*N_U + 3*N_S + 3*N_D + N_P;
  
  N_Active_U = SqmatrixA11->GetActiveBound();
  N_Active_S = SqmatrixG11->GetActiveBound();
  N_Active_D = SqmatrixH11->GetActiveBound();

  KColA = SqmatrixA11->GetKCol();
  RowPtrA = SqmatrixA11->GetRowPtr();
  
  KColG = SqmatrixG11->GetKCol();
  RowPtrG = SqmatrixG11->GetRowPtr();
  
  KColH = SqmatrixH11->GetKCol();
  RowPtrH = SqmatrixH11->GetRowPtr();

  KColB = MatrixB1->GetKCol();
  RowPtrB = MatrixB1->GetRowPtr();

  KColBT = MatrixB1T->GetKCol();
  RowPtrBT = MatrixB1T->GetRowPtr();
  
  KColC = MatrixC11->GetKCol();
  RowPtrC = MatrixC11->GetRowPtr();
  
  KColD = MatrixD11->GetKCol();
  RowPtrD = MatrixD11->GetRowPtr();
  
  KColE = MatrixE11->GetKCol();
  RowPtrE = MatrixE11->GetRowPtr();
  
  KColJ = MatrixJ11->GetKCol();
  RowPtrJ = MatrixJ11->GetRowPtr();

  EntriesA11 = SqmatrixA11->GetEntries();
  EntriesA12 = SqmatrixA12->GetEntries();
  EntriesA21 = SqmatrixA21->GetEntries();
  EntriesA22 = SqmatrixA22->GetEntries();
  
  EntriesG11 = SqmatrixG11->GetEntries();
  EntriesG12 = SqmatrixG12->GetEntries();
  EntriesG21 = SqmatrixG21->GetEntries();
  EntriesG22 = SqmatrixG22->GetEntries();
  EntriesG23 = SqmatrixG23->GetEntries();
  EntriesG32 = SqmatrixG32->GetEntries();
  EntriesG33 = SqmatrixG33->GetEntries();
  
  EntriesH11 = SqmatrixH11->GetEntries();
  EntriesH22 = SqmatrixH22->GetEntries();
  EntriesH33 = SqmatrixH33->GetEntries();
  
  EntriesB1 = MatrixB1->GetEntries();
  EntriesB2 = MatrixB2->GetEntries();
  EntriesB1T = MatrixB1T->GetEntries();
  EntriesB2T = MatrixB2T->GetEntries();
  
  EntriesC11 = MatrixC11->GetEntries();
  EntriesC12 = MatrixC12->GetEntries();
  EntriesC22 = MatrixC22->GetEntries();
  EntriesC23 = MatrixC23->GetEntries();
  
  EntriesD11 = MatrixD11->GetEntries();
  EntriesD12 = MatrixD12->GetEntries();
  EntriesD21 = MatrixD21->GetEntries();
  EntriesD22 = MatrixD22->GetEntries();
  EntriesD31 = MatrixD31->GetEntries();
  EntriesD32 = MatrixD32->GetEntries();
  
  EntriesE11 = MatrixE11->GetEntries();
  EntriesE12 = MatrixE12->GetEntries();
  EntriesE22 = MatrixE22->GetEntries();
  EntriesE23 = MatrixE23->GetEntries();
  
  EntriesJ11 = MatrixJ11->GetEntries();
  EntriesJ21 = MatrixJ21->GetEntries();
  EntriesJ22 = MatrixJ22->GetEntries();
  EntriesJ32 = MatrixJ32->GetEntries();

  N_Entries = 4*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U] + 7*RowPtrG[N_S] + 3*RowPtrH[N_D] + 4*RowPtrC[N_U] + 4*RowPtrE[N_U] + 6*RowPtrD[N_S] + 4*RowPtrJ[N_D];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_Tot+1];
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

      Entries[pos] = (i<N_Active_U)?EntriesA12[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active_U)
    {
      begin = RowPtrC[i];
      end = RowPtrC[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesC11[j];
        KCol[pos] = KColC[j]+2*N_U;
        pos++;
	
	Entries[pos] = EntriesC12[j];
        KCol[pos] = KColC[j]+2*N_U + N_S;
        pos++;
      }         
            
      begin = RowPtrE[i];
      end = RowPtrE[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesE11[j];
        KCol[pos] = KColE[j]+2*N_U + 3*N_S;
        pos++;
	
	Entries[pos] = EntriesE12[j];
        KCol[pos] = KColE[j]+2*N_U + 3*N_S + N_D;
        pos++;
      }
     
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1T[j];
        KCol[pos] = KColBT[j]+2*N_U + 3*N_S + 3*N_D;
        pos++;
      }
      
    } // if(i<N_Active_U)
    
    RowPtr[i+1] = pos;
  } // for(i=0;i<N_U;i++)
 
  
  for(i=0;i<N_U;i++)
  { 
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active_U)?EntriesA21[j]:0;
      KCol[pos] = KColA[j];
      pos++;

      Entries[pos] = EntriesA22[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active_U)
    {
     begin = RowPtrC[i];
      end = RowPtrC[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesC22[j];
        KCol[pos] = KColC[j]+2*N_U + N_S;
        pos++;
	
	Entries[pos] = EntriesC23[j];
        KCol[pos] = KColC[j]+2*N_U + 2*N_S;
        pos++;
      }       
            
      begin = RowPtrE[i];
      end = RowPtrE[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesE22[j];
        KCol[pos] = KColE[j]+2*N_U + 3*N_S + N_D;
        pos++;
	
	Entries[pos] = EntriesE23[j];
        KCol[pos] = KColE[j]+2*N_U + 3*N_S + 2*N_D;
        pos++;
      }
      
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2T[j];
        KCol[pos] = KColBT[j]+2*N_U + 3*N_S + 3*N_D;
        pos++;
      }
      
    } // if(i<N_Active_U)
    
    RowPtr[N_U+i+1] = pos;
  }

    for(i=0;i<N_S;i++)
  {
    if(i<N_Active_S)
    {
      begin = RowPtrD[i];
      end = RowPtrD[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesD11[j];
        KCol[pos] = KColD[j];
        pos++;
	
	Entries[pos] = EntriesD12[j];
        KCol[pos] = KColD[j]+ N_U;
        pos++;
      }         
      
    } // if(i<N_Active_S)
    
    
    begin = RowPtrG[i];
    end = RowPtrG[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesG11[j];
      KCol[pos] = KColG[j] + 2*N_U;
      pos++;
      Entries[pos] = (i<N_Active_S)?EntriesG12[j]:0;
      KCol[pos] = KColG[j] + 2*N_U + N_S;
      pos++;
    }
  
    RowPtr[2*N_U+i+1] = pos;
  }
  
      for(i=0;i<N_S;i++)
  {
    if(i<N_Active_S)
    {
      begin = RowPtrD[i];
      end = RowPtrD[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesD21[j];
        KCol[pos] = KColD[j];
        pos++;
	
	Entries[pos] = EntriesD22[j];
        KCol[pos] = KColD[j]+ N_U;
        pos++;
      }         
      
    } // if(i<N_Active_S)
    
    
    begin = RowPtrG[i];
    end = RowPtrG[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active_S)?EntriesG21[j]:0;
      KCol[pos] = KColG[j] + 2*N_U;
      pos++;
      Entries[pos] = EntriesG22[j];
      KCol[pos] = KColG[j] + 2*N_U + N_S;
      pos++;
      Entries[pos] = (i<N_Active_S)?EntriesG23[j]:0;
      KCol[pos] = KColG[j] + 2*N_U + 2*N_S;
      pos++;
    }
  
    RowPtr[2*N_U+N_S+i+1] = pos;
  }
  
    for(i=0;i<N_S;i++)
  {
    if(i<N_Active_S)
    {
      begin = RowPtrD[i];
      end = RowPtrD[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesD31[j];
        KCol[pos] = KColD[j];
        pos++;
	
	Entries[pos] = EntriesD32[j];
        KCol[pos] = KColD[j]+ N_U;
        pos++;
      }         
      
    } // if(i<N_Active_S)
    
    
    begin = RowPtrG[i];
    end = RowPtrG[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active_S)?EntriesG32[j]:0;
      KCol[pos] = KColG[j] + 2*N_U + N_S;
      pos++;
      Entries[pos] = EntriesG33[j];
      KCol[pos] = KColG[j] + 2*N_U + 2*N_S;
      pos++;
    }
  
    RowPtr[2*N_U+2*N_S+i+1] = pos;
  }
  
   for(i=0;i<N_D;i++)
  {
    if(i<N_Active_D)
    {
      begin = RowPtrJ[i];
      end = RowPtrJ[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesJ11[j];
        KCol[pos] = KColJ[j];
        pos++;
      }         
      
    } // if(i<N_Active_D)
      
    begin = RowPtrH[i];
    end = RowPtrH[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesH11[j];
      KCol[pos] = KColH[j] + 2*N_U + 3*N_S;
      pos++;
    }
  
    RowPtr[2*N_U+3*N_S+i+1] = pos;
  }
 
   for(i=0;i<N_D;i++)
  {
    if(i<N_Active_D)
    {
      begin = RowPtrJ[i];
      end = RowPtrJ[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesJ21[j];
        KCol[pos] = KColJ[j];
        pos++;
	Entries[pos] = EntriesJ22[j];
        KCol[pos] = KColJ[j] + N_U;
        pos++;
      }         
      
    } // if(i<N_Active_D)
      
    begin = RowPtrH[i];
    end = RowPtrH[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesH22[j];
      KCol[pos] = KColH[j] + 2*N_U + 3*N_S + N_D;
      pos++;
    }
  
    RowPtr[2*N_U+3*N_S+N_D+i+1] = pos;
  }
 
   for(i=0;i<N_D;i++)
  {
    if(i<N_Active_D)
    {
      begin = RowPtrJ[i];
      end = RowPtrJ[i+1];
      for(j=begin;j<end;j++)
      {
	Entries[pos] = EntriesJ32[j];
        KCol[pos] = KColJ[j] + N_U;
        pos++;
      }         
      
    } // if(i<N_Active_D)
      
    begin = RowPtrH[i];
    end = RowPtrH[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesH33[j];
      KCol[pos] = KColH[j] + 2*N_U + 3*N_S + 2*N_D;
      pos++;
    }
  
    RowPtr[2*N_U+3*N_S+2*N_D+i+1] = pos;
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
    RowPtr[2*N_U+3*N_S+3*N_D+i+1] = pos;
  }

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[2*N_U+3*N_S+3*N_D];
    end = RowPtr[2*N_U+3*N_S+3*N_D+1];
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 2*N_U+3*N_S+3*N_D;
    rhs[2*N_U+3*N_S+3*N_D] = 0;
  }

  // sort matrix
  for(i=0;i<N_Tot;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  
  t2 = GetTime();

  ret = umfpack_di_symbolic(N_Tot, N_Tot, RowPtr, KCol, Entries, &Symbolic, null, null);
  if(ret!=0)
  OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)  
  OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)    
  OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }
  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}



// NSE(Type-4)-CST-2D (Monolithic) LPS
void DirectSolver(TSquareMatrix2D *SqmatrixA11, TSquareMatrix2D *SqmatrixA12, TSquareMatrix2D *SqmatrixA21, TSquareMatrix2D *SqmatrixA22,
		  TSquareMatrix2D *SqmatrixG11, TSquareMatrix2D *SqmatrixG12, TSquareMatrix2D *SqmatrixG21 , TSquareMatrix2D *SqmatrixG22, TSquareMatrix2D *SqmatrixG23 , TSquareMatrix2D *SqmatrixG32 , TSquareMatrix2D *SqmatrixG33,
		  TMatrix2D *MatrixB1,  TMatrix2D *MatrixB2, TMatrix2D *MatrixB1T, TMatrix2D *MatrixB2T, 
                  TMatrix2D *MatrixC11, TMatrix2D *MatrixC12, TMatrix2D *MatrixC22, TMatrix2D *MatrixC23,
		  TMatrix2D *MatrixD11, TMatrix2D *MatrixD12, TMatrix2D *MatrixD21, TMatrix2D *MatrixD22, TMatrix2D *MatrixD31, TMatrix2D *MatrixD32,
                  double *rhs, double *sol)
{
  int *KColA, *RowPtrA, *KColB, *RowPtrB, *KColBT, *RowPtrBT, *KColG, *RowPtrG;
  int *KColC, *RowPtrC, *KColD, *RowPtrD;
  
  double *EntriesA11, *EntriesA12, *EntriesA21, *EntriesA22;
  double *EntriesG11, *EntriesG12, *EntriesG21, *EntriesG22, *EntriesG23, *EntriesG32, *EntriesG33;
  double *EntriesB1, *EntriesB2, *EntriesB1T, *EntriesB2T;
  double *EntriesC11, *EntriesC12, *EntriesC22, *EntriesC23;
  double *EntriesD11, *EntriesD12, *EntriesD21, *EntriesD22, *EntriesD31, *EntriesD32;
  
  int N_Tot, N_U, N_P, N_S, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active_U, N_Active_S;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = SqmatrixA11->GetN_Rows();
  N_S = SqmatrixG11->GetN_Rows();
  N_P = MatrixB1->GetN_Rows();
  N_Tot = 2*N_U + 3*N_S + N_P;
  
  N_Active_U = SqmatrixA11->GetActiveBound();
  N_Active_S = SqmatrixG11->GetActiveBound();

  KColA = SqmatrixA11->GetKCol();
  RowPtrA = SqmatrixA11->GetRowPtr();
  
  KColG = SqmatrixG11->GetKCol();
  RowPtrG = SqmatrixG11->GetRowPtr();

  KColB = MatrixB1->GetKCol();
  RowPtrB = MatrixB1->GetRowPtr();

  KColBT = MatrixB1T->GetKCol();
  RowPtrBT = MatrixB1T->GetRowPtr();
  
  KColC = MatrixC11->GetKCol();
  RowPtrC = MatrixC11->GetRowPtr();
  
  KColD = MatrixD11->GetKCol();
  RowPtrD = MatrixD11->GetRowPtr();

  EntriesA11 = SqmatrixA11->GetEntries();
  EntriesA12 = SqmatrixA12->GetEntries();
  EntriesA21 = SqmatrixA21->GetEntries();
  EntriesA22 = SqmatrixA22->GetEntries();
  
  EntriesG11 = SqmatrixG11->GetEntries();
  EntriesG12 = SqmatrixG12->GetEntries();
  EntriesG21 = SqmatrixG21->GetEntries();
  EntriesG22 = SqmatrixG22->GetEntries();
  EntriesG23 = SqmatrixG23->GetEntries();
  EntriesG32 = SqmatrixG32->GetEntries();
  EntriesG33 = SqmatrixG33->GetEntries();
    
  EntriesB1 = MatrixB1->GetEntries();
  EntriesB2 = MatrixB2->GetEntries();
  EntriesB1T = MatrixB1T->GetEntries();
  EntriesB2T = MatrixB2T->GetEntries();
  
  EntriesC11 = MatrixC11->GetEntries();
  EntriesC12 = MatrixC12->GetEntries();
  EntriesC22 = MatrixC22->GetEntries();
  EntriesC23 = MatrixC23->GetEntries();
  
  EntriesD11 = MatrixD11->GetEntries();
  EntriesD12 = MatrixD12->GetEntries();
  EntriesD21 = MatrixD21->GetEntries();
  EntriesD22 = MatrixD22->GetEntries();
  EntriesD31 = MatrixD31->GetEntries();
  EntriesD32 = MatrixD32->GetEntries();

  N_Entries = 4*RowPtrA[N_U] + 2*RowPtrB[N_P] + 2*RowPtrBT[N_U] + 7*RowPtrG[N_S] + 4*RowPtrC[N_U] + 6*RowPtrD[N_S];
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_Tot+1];
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

      Entries[pos] = (i<N_Active_U)?EntriesA12[j]:0;
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active_U)
    {
      begin = RowPtrC[i];
      end = RowPtrC[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesC11[j];
        KCol[pos] = KColC[j]+2*N_U;
        pos++;
	
	Entries[pos] = EntriesC12[j];
        KCol[pos] = KColC[j]+2*N_U + N_S;
        pos++;
      }         
                 
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB1T[j];
        KCol[pos] = KColBT[j]+2*N_U + 3*N_S;
        pos++;
      }
      
    } // if(i<N_Active_U)
    
    RowPtr[i+1] = pos;
  } // for(i=0;i<N_U;i++)
 
  
  for(i=0;i<N_U;i++)
  {
    begin = RowPtrA[i];
    end = RowPtrA[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active_U)?EntriesA21[j]:0;
      KCol[pos] = KColA[j];
      pos++;

      Entries[pos] = EntriesA22[j];
      KCol[pos] = KColA[j]+N_U;
      pos++;
    }

    if(i<N_Active_U)
    {
      begin = RowPtrC[i];
      end = RowPtrC[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesC22[j];
        KCol[pos] = KColC[j]+2*N_U + N_S;
        pos++;
	
	Entries[pos] = EntriesC23[j];
        KCol[pos] = KColC[j]+2*N_U + 2*N_S;
        pos++;
      }         
                 
      begin = RowPtrBT[i];
      end = RowPtrBT[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesB2T[j];
        KCol[pos] = KColBT[j]+2*N_U + 3*N_S;
        pos++;
      }
      
    } // if(i<N_Active_U)
    
    RowPtr[N_U+i+1] = pos;
  }

    for(i=0;i<N_S;i++)
  {
    if(i<N_Active_S)
    {
      begin = RowPtrD[i];
      end = RowPtrD[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesD11[j];
        KCol[pos] = KColD[j];
        pos++;
	
	Entries[pos] = EntriesD12[j];
        KCol[pos] = KColD[j]+ N_U;
        pos++;
      }         
      
    } // if(i<N_Active_S)
    
    
    begin = RowPtrG[i];
    end = RowPtrG[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = EntriesG11[j];
      KCol[pos] = KColG[j] + 2*N_U;
      pos++;
      Entries[pos] = (i<N_Active_S)?EntriesG12[j]:0;
      KCol[pos] = KColG[j] + 2*N_U + N_S;
      pos++;
    }
  
    RowPtr[2*N_U+i+1] = pos;
  }
  
      for(i=0;i<N_S;i++)
  {
    if(i<N_Active_S)
    {
      begin = RowPtrD[i];
      end = RowPtrD[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesD21[j];
        KCol[pos] = KColD[j];
        pos++;
	
	Entries[pos] = EntriesD22[j];
        KCol[pos] = KColD[j]+ N_U;
        pos++;
      }         
      
    } // if(i<N_Active_S)
    
    
    begin = RowPtrG[i];
    end = RowPtrG[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active_S)?EntriesG21[j]:0;
      KCol[pos] = KColG[j] + 2*N_U;
      pos++;
      Entries[pos] = EntriesG22[j];
      KCol[pos] = KColG[j] + 2*N_U + N_S;
      pos++;
      Entries[pos] = (i<N_Active_S)?EntriesG23[j]:0;
      KCol[pos] = KColG[j] + 2*N_U + 2*N_S;
      pos++;
    }
  
    RowPtr[2*N_U+N_S+i+1] = pos;
  }
  
    for(i=0;i<N_S;i++)
  {
    if(i<N_Active_S)
    {
      begin = RowPtrD[i];
      end = RowPtrD[i+1];
      for(j=begin;j<end;j++)
      {
        Entries[pos] = EntriesD31[j];
        KCol[pos] = KColD[j];
        pos++;
	
	Entries[pos] = EntriesD32[j];
        KCol[pos] = KColD[j]+ N_U;
        pos++;
      }         
      
    } // if(i<N_Active_S)
      
    begin = RowPtrG[i];
    end = RowPtrG[i+1];
    for(j=begin;j<end;j++)
    {
      Entries[pos] = (i<N_Active_S)?EntriesG32[j]:0;
      KCol[pos] = KColG[j] + 2*N_U + N_S;
      pos++;
      Entries[pos] = EntriesG33[j];
      KCol[pos] = KColG[j] + 2*N_U + 2*N_S;
      pos++;
    }
  
    RowPtr[2*N_U+2*N_S+i+1] = pos;
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
    RowPtr[2*N_U+3*N_S+i+1] = pos;
  }

  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[2*N_U+3*N_S];
    end = RowPtr[2*N_U+3*N_S+1];
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 2*N_U+3*N_S;
    rhs[2*N_U+3*N_S] = 0;
  }

  // sort matrix
  for(i=0;i<N_Tot;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_Tot;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */
  
//   double sum1 = 0.0, sum2 = 0.0;
//     for(i=0;i<N_P;i++)
//   {
//     begin = RowPtrB[i];
//     end = RowPtrB[i+1];
//     for(j=begin;j<end;j++)
//     {
//       sum1 += EntriesB1[j];
//       sum2 += EntriesB2[j];
//     }
//   }
//   cout<<"B1 : "<<sum1<<"\n"<<"B2 : "<<sum2<<"\n";
//   
//   sum1= 0.0; 
//   sum2 = 0.0;
//   
//       for(i=0;i<N_U;i++)
//   {
//     begin = RowPtrBT[i];
//     end = RowPtrBT[i+1];
//     for(j=begin;j<end;j++)
//     {
//       sum1 += EntriesB1T[j];
//       sum2 += EntriesB2T[j];
//     }
//   }
//   cout<<"B1T : "<<sum1<<"\n"<<"B2T : "<<sum2<<"\n";
//   
//    double sum11= 0.0, sum12=0.0, sum21=0.0, sum22=0.0; 
//   
//         for(i=0;i<N_U;i++)
//   {
//     begin = RowPtrA[i];
//     end = RowPtrA[i+1];
//     for(j=begin;j<end;j++)
//     {
//       sum11 += EntriesA11[j];
//       sum12 += EntriesA12[j];
//       sum21 += EntriesA21[j];
//       sum22 += EntriesA22[j];
//     }
//   }
//   cout<<"A11 : "<<sum11<<"\nA12 : "<<sum12<<"\nA21 : "<<sum21<<"\nA22 : "<<sum22<<"\n";
//   
//   sum11= 0.0;
//   sum12=0.0;
//   sum21=0.0;
//   sum22=0.0; 
//     for(i=0;i<N_U;i++)
//   {
//     begin = RowPtrC[i];
//     end = RowPtrC[i+1];
//     for(j=begin;j<end;j++)
//     {
//       sum11 += EntriesC11[j];
//       sum12 += EntriesC12[j];
//       sum21 += EntriesC22[j];
//       sum22 += EntriesC23[j];
//     }
//   }
//   cout<<"C11 : "<<sum11<<"\nC12 : "<<sum12<<"\nC22 : "<<sum21<<"\nC23 : "<<sum22<<"\n";
//   
//   
//   
//   
//    cout<<"sol u: "<<Ddot(2*N_U,sol, sol)<<"\n";
// //  cout<<"sol s: "<<Ddot(3*N_S,sol + 2*N_U, sol+2*N_U)<<"\n";
//  cout<<"sol p: "<<Ddot(N_P,sol+2*N_U+3*N_S, sol+2*N_U+3*N_S)<<"\n";
//   
//  cout<<"rhs u: "<<Ddot(2*N_U,rhs, rhs)<<"\n";
// //  cout<<"rhs s: "<<Ddot(3*N_S,rhs + 2*N_U, rhs+2*N_U)<<"\n";
//  cout<<"rhs p: "<<Ddot(N_P,rhs+2*N_U+3*N_S, rhs+2*N_U+3*N_S)<<"\n";
 
  t2 = GetTime();

  ret = umfpack_di_symbolic(N_Tot, N_Tot, RowPtr, KCol, Entries, &Symbolic, null, null);
  if(ret!=0)
  OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)  
  OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)    
  OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();
  
//  cout<<"sol u: "<<Ddot(2*N_U,sol, sol)<<"\n";
// //  cout<<"sol s: "<<Ddot(3*N_S,sol + 2*N_U, sol+2*N_U)<<"\n";
//  cout<<"sol p: "<<Ddot(N_P,sol+2*N_U+3*N_S, sol+2*N_U+3*N_S)<<"\n";
//   
//  cout<<"rhs u: "<<Ddot(2*N_U,rhs, rhs)<<"\n";
// //  cout<<"rhs s: "<<Ddot(3*N_S,rhs + 2*N_U, rhs+2*N_U)<<"\n";
//  cout<<"rhs p: "<<Ddot(N_P,rhs+2*N_U+3*N_S, rhs+2*N_U+3*N_S)<<"\n";
//  exit(0);
  
  
  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

  if(TDatabase::ParamDB->SC_VERBOSE>1)
  {
     cout << "UMFPACK Time:";
     cout << "  data prep: " << t2-t1 << "s ";
     cout << "  symbolic: " << t3-t2 << "s ";
     cout << "  numeric: " << t4-t3 << "s ";
     cout << "  solve: " << t5-t4 << "s "<< endl;
     cout << "UMFPACK total time: " << t5-t1 << "s "<< endl;
  }
  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}






#ifdef __3D__
//****************************************************************************/
//
// for NSTYPE == 2
//
//****************************************************************************/

void DirectSolver(TSquareMatrix3D *sqmatrixA,
TMatrix3D *matrixB1T, TMatrix3D *matrixB2T,
TMatrix3D *matrixB3T,
TMatrix3D *matrixB1,  TMatrix3D *matrixB2,
TMatrix3D *matrixB3,
double *rhs, double *sol)
{
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA;
  double *EntriesB1, *EntriesB2, *EntriesB3;
  double *EntriesB1T, *EntriesB2T, *EntriesB3T;
  int N_, N_U, N_P, N_Entries;
  double *Entries;
  int *KCol, *RowPtr;
  double *null = (double *) NULL;
  void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;

  t1 = GetTime();
  N_U = sqmatrixA->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 3*N_U + N_P;
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
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
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

double len = 0;
  if(TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE)
  {
    // pressure constant
    begin = RowPtr[3*N_U];
    end = RowPtr[3*N_U+1];
    
    len=end- begin;
//     cout << "No. Internal Pressure row entries " <<len-- <<endl;
    
    for(j=begin+1;j<end;j++)
      Entries[j] = 0;
    Entries[begin] = 1;
    KCol[begin] = 3*N_U;
    rhs[3*N_U] = 0;
  }

//   cout << "Total entries : " <<  RowPtr[N_]-len << endl;
  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i


//   for(i=0;i<N_;i++)
//   {
//     cout << "=====" << RowPtr[i] << " " << RowPtr[i+1] << endl;
//     for(j=RowPtr[i];j<RowPtr[i+1];j++)
//       cout << i << " " << KCol[j] << " " << Entries[j] << endl;
//   }

//   for(i=0;i<N_U;i++)
//   cout << "Rhs " << i << " " << rhs[i] <<" " << rhs[i+1] <<  " " << rhs[i+2] << endl;

//   for(i=0;i<N_P;i++)
//   cout << "Rhs " << i << " " << rhs[3*N_U+i] << endl;
// 

// exit(0);
  

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  if(ret!=0)
  OutPut("symbolic: " << ret << endl);
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if(ret!=0)
   OutPut("numeric: " << ret << endl);
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if(ret!=0)
   OutPut("solve: " << ret << endl);
  umfpack_di_free_numeric(&Numeric);
  t5 = GetTime();

  delete [] Entries;
  delete [] KCol;
  delete [] RowPtr;

//   cout << "UMFPACK:";
//   cout << "  data prep: " << t2-t1 << " ";
//   cout << "  symbolic: " << t3-t2 << " ";
//   cout << "  numeric: " << t4-t3 << " ";
//   cout << "  solve: " << t5-t4 << endl;
//   cout << "UMFPACK total time: " << t5-t1 << endl;
}


//****************************************************************************/
//
// for NSTYPE == 4
// flag = 0 ==> allocate and LU factorization, solve
// flag = 1 ==> only forward and backward solve
// flag = 2 ==> free memory
// flag = 3 ==> allocate, LU factorization, solve  and free memory
// flag = 4 ==> free memory
//
//****************************************************************************/

void DirectSolver(TSquareMatrix3D *sqmatrixA11, TSquareMatrix3D *sqmatrixA12,
TSquareMatrix3D *sqmatrixA13,
TSquareMatrix3D *sqmatrixA21, TSquareMatrix3D *sqmatrixA22,
TSquareMatrix3D *sqmatrixA23,
TSquareMatrix3D *sqmatrixA31, TSquareMatrix3D *sqmatrixA32,
TSquareMatrix3D *sqmatrixA33,
TMatrix3D *matrixB1T, TMatrix3D *matrixB2T, TMatrix3D *matrixB3T,
TMatrix3D *matrixB1,  TMatrix3D *matrixB2, TMatrix3D *matrixB3,
double *rhs, double *sol, int flag)
{
  if (TDatabase::ParamDB->SC_VERBOSE>3)  
      OutPut("umf3d"<<endl);
  int *KColA, *RowPtrA;
  int *KColB, *RowPtrB;
  int *KColBT, *RowPtrBT;
  double *EntriesA11, *EntriesA12, *EntriesA13, *EntriesA21;
  double *EntriesA22, *EntriesA23, *EntriesA31, *EntriesA32, *EntriesA33;
  double *EntriesB1, *EntriesB2,  *EntriesB3, *EntriesB1T, *EntriesB2T, *EntriesB3T;
  int N_, N_U, N_P, N_Entries;
  static double *Entries;
  static int *KCol, *RowPtr;
  double *null = (double *) NULL;
  static void *Symbolic, *Numeric;
  int i, j, k, l, begin, end, ret, pos;
  double value;
  int N_Active;
  double t1, t2, t3, t4, t5;
  int verbose = TDatabase::ParamDB->SC_VERBOSE;

  if (flag==4)
  {
    umfpack_di_free_numeric(&Numeric);

    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
    return;
  }

  if (verbose>3)
  {
      OutPut("flag: " << flag << endl);
  }

  t1 = GetTime();
  if (flag==0 || flag==3)
  {
  N_U = sqmatrixA11->GetN_Rows();
  N_P = matrixB1->GetN_Rows();
  N_ = 3*N_U + N_P;
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
  Entries = new double[N_Entries];
  KCol = new int[N_Entries];

  RowPtr = new int[N_+1];
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
    rhs[3*N_U] = 0;
  }
  
  // sort matrix
  for(i=0;i<N_;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(j=begin;j<end;j++)
    {
      for(k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }                              // endfor i

  /*
  for(i=0;i<N_;i++)
  {
    for(j=RowPtr[i];j<RowPtr[i+1];j++)
      cout << i << " " << KCol[j] << " " << Entries[j] << endl;
  }
  */

  t2 = GetTime();

  ret = umfpack_di_symbolic(N_, N_, RowPtr, KCol, Entries, &Symbolic, null, null);
  if (ret!=0)
  {
      OutPut("WARNING: symbolic: " << ret << endl);
  }
  t3 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, null, null);
  if (ret!=0)
  {
      OutPut("WARNING: numeric: " << ret << endl);
  }
  t4 = GetTime();
  umfpack_di_free_symbolic(&Symbolic);
  }
  t4 = GetTime();
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries,
    sol, rhs, Numeric, null, null);
  if (ret!=0)
  {
      OutPut("WARNING: solve: " << ret << endl);
  }
  t5 = GetTime();

  if (flag==2 || flag==3)
  {
      umfpack_di_free_numeric(&Numeric);
      delete [] Entries;
      delete [] KCol;
      delete [] RowPtr;
  }

  if (verbose>3)
  {
  OutPut( "UMFPACK:");
  if (flag==0 || flag==3)
  {
      OutPut( "  data prep: " << t2-t1 << " ");
      OutPut( "  symbolic: " << t3-t2 << " ");
      OutPut( "  numeric: " << t4-t3 << " ");
  }
  OutPut( "  solve: " << t5-t4);
  OutPut( "  total time: " << t5-t1 << endl);
  }
  /*
  for(i=0;i<N_;i++)
    cout << setw(6) << i << setw(30) << sol[i] << endl;
  */
}



void DirectSolver(TSquareMatrix3D **sqmatrices, int n_row, int n_column,
                  double *sol, double *rhs)
{
  int N_Active;
  int N_Rows;
  int N_Entries, N_Row, begin, end, pos, l;
  int *RowPtr, *rowptr;
  int *KCol, *kcol;
  double *Entries, *entries, value;
  double t1, t2;
  
  N_Active = sqmatrices[0]->GetActiveBound();
  N_Rows   = sqmatrices[0]->GetN_Rows();
  
  N_Entries = n_row*n_column*sqmatrices[0]->GetN_Entries();
  
  Entries = new double [N_Entries];
  RowPtr  = new int [N_Rows*n_row+1];
  KCol    = new int [N_Entries];
  
  N_Row = N_Rows*n_row;
  
  pos = 0;
  RowPtr[0] = 0;
  
  for (int i=0;i<n_row;++i)
  {
    for (int row=0;row<N_Rows;++row)
    {
      for (int j=0;j<n_column;++j)
      {
        if ( i != j && row >= N_Active ) continue;
 
        rowptr = sqmatrices[i*n_column+j]->GetRowPtr();
        kcol   = sqmatrices[i*n_column+j]->GetKCol();
        entries = sqmatrices[i*n_column+j]->GetEntries();
 
        begin = rowptr[row];
        end   = rowptr[row+1];

        for (int loc_pos=begin;loc_pos<end;++loc_pos)
         {
          Entries[pos] = entries[loc_pos];
          KCol[pos] = kcol[loc_pos] + j*N_Rows;
          ++pos;
         }
      }
      RowPtr[i*N_Rows+row+1] = pos;
    }
  }
  
   // sort matrix
  for(int i=0;i<N_Row;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(int j=begin;j<end;j++)
    {
      for(int k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }       
  
  void *Symbolic, *Numeric;
  int ret;
  
  ret = umfpack_di_symbolic(N_Row, N_Row, RowPtr, KCol, Entries, &Symbolic, NULL, NULL);
  if ( ret != 0) 
  {
    OutPut("symbolic: " << ret << endl);
    exit(0);
  }
//   t1 = GetTime();
  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, NULL, NULL);
  if ( ret != 0 )
  {
    OutPut("numeric: " << ret << endl);
    exit(0);
  }
//   t2 = GetTime();
//   OutPut("numeric: " << ret << " "  << t2-t1 << endl);
    
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries, sol, rhs, Numeric, NULL, NULL);
  if ( ret != 0) 
  {
    OutPut("solve: " << ret << endl);
    exit(0);
  }
//   t1 = GetTime();
//   OutPut("numeric: " << ret << " "  << t1-t2 << endl);  
  
  umfpack_di_free_symbolic(&Symbolic);
  umfpack_di_free_numeric(&Numeric);
  
  delete [] Entries;
  delete [] RowPtr;
  delete [] KCol;
  
}

// rb_flag = 0 ==> allocation and LU-decomposition forward/backward.
// rb_flag = 1 ==> only forward/backward.
// rb_flag = 2 ==> forward/backward and free up memory
// rb_flag = 3 ==> allocation, LU-decomposition, forward/backward, free up memory
// rb_flag = 4 ==> only free up memory
void DirectSolver(TSquareMatrix3D **sqmatrices, int n_row, int n_column,
                  double *sol, double *rhs, double *&Entries,
                   int *&KCol, int *&RowPtr, void *&Symbolic, void *&Numeric, int rb_flag)
{
  int N_Active, ret;
  int N_Rows;
  int N_Entries, N_Row, begin, end, pos, l;
  int *rowptr;
  int  *kcol;
  double *entries, value;
  double t1, t2;
  
 if (rb_flag==4)
  {
   umfpack_di_free_numeric(&Numeric);
 
    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
    return;
  }    
 
  N_Active = sqmatrices[0]->GetActiveBound();
  N_Rows   = sqmatrices[0]->GetN_Rows();
  N_Entries = n_row*n_column*sqmatrices[0]->GetN_Entries();
  
 if (TDatabase::ParamDB->SC_VERBOSE>=3)
  {
   OutPut("rb_flag: " << rb_flag << endl);
  }
  

 if (rb_flag==0 || rb_flag==3)
 {   
  Entries = new double [N_Entries];
  RowPtr  = new int [N_Rows*n_row+1];
  KCol    = new int [N_Entries];
  
  N_Row = N_Rows*n_row;
  
  pos = 0;
  RowPtr[0] = 0;
  
  for (int i=0;i<n_row;++i)
   {
    for (int row=0;row<N_Rows;++row)
    {
      for (int j=0;j<n_column;++j)
      {
        if ( i != j && row >= N_Active ) continue;
 
        rowptr = sqmatrices[i*n_column+j]->GetRowPtr();
        kcol   = sqmatrices[i*n_column+j]->GetKCol();
        entries = sqmatrices[i*n_column+j]->GetEntries();
 
        begin = rowptr[row];
        end   = rowptr[row+1];

        for (int loc_pos=begin;loc_pos<end;++loc_pos)
         {
          Entries[pos] = entries[loc_pos];
          KCol[pos] = kcol[loc_pos] + j*N_Rows;
          ++pos;
         }
      }
      RowPtr[i*N_Rows+row+1] = pos;
    }
  }
  
   // sort matrix
  for(int i=0;i<N_Row;i++)
  {
    begin=RowPtr[i];
    end=RowPtr[i+1];

    for(int j=begin;j<end;j++)
    {
      for(int k=j+1;k<end;k++)
      {
        if(KCol[j] > KCol[k])
        {
          l = KCol[j];      value = Entries[j];
          KCol[j] = KCol[k]; Entries[j] = Entries[k];
          KCol[k] = l;       Entries[k] = value;
        }                        // endif
      }                          // endfor k
    }                            // endfor j
  }       

  ret = umfpack_di_symbolic(N_Row, N_Row, RowPtr, KCol, Entries, &Symbolic, NULL, NULL);
  UMFPACK_return(ret);

  ret = umfpack_di_numeric(RowPtr, KCol, Entries, Symbolic, &Numeric, NULL, NULL);
  UMFPACK_return(ret);  
  umfpack_di_free_symbolic(&Symbolic);
 } //if (rb_flag==0 || rb_flag==3)
 
  ret = umfpack_di_solve(UMFPACK_At, RowPtr, KCol, Entries, sol, rhs, Numeric, NULL, NULL);
  UMFPACK_return(ret);
 
 if (rb_flag==2 || rb_flag==3)
  {
   umfpack_di_free_numeric(&Numeric);
    
    delete [] Entries;
    delete [] KCol;
    delete [] RowPtr;
  }  
}

#endif
